///
/// \file ThermSim.hh
///


#pragma once

#include <TLorentzVector.h>
#include <TH2D.h>

#include <vector>
#include <thread>


struct Particle {
  TLorentzVector x;
  TLorentzVector p;
  int charge;

  Particle(TLorentzVector X, TLorentzVector P)
    : x(X)
    , p(P)
    { }

  double pt() const
    { return p.Pt(); }
};


struct Event {
  using container_t = std::vector<Particle>;

  container_t particles;

  Event()
    { }

  Event(Event &&e)
    : particles(std::move(e.particles))
    { }

  Event(container_t &&p)
    : particles(std::move(p))
    { }

  Event(const container_t &p)
    : particles(p)
    { }

  Event& operator=(Event &&ev)
    {
      particles = std::move(ev.particles);
      return *this;
    }

  container_t::iterator begin()
    { return particles.begin(); }

  container_t::iterator end()
    { return particles.end(); }

};

struct FemtoRunner {
  std::mutex &mutex;
  std::deque<std::unique_ptr<Event>> &event_queue;
  std::condition_variable &cv;

  std::deque<Event> mixing_queue;

  std::unique_ptr<TH2D>
    num,
    den;

  FemtoRunner(std::deque<std::unique_ptr<Event>> &q,
              std::mutex &m,
              std::condition_variable &c)
    : mutex(m)
    , event_queue(q)
    , cv(c)
    , mixing_queue()
    , num(std::make_unique<TH2D>("num",
                                 "Numerator; q_{inv}; k_{T} (GeV);",
                                 200, 0.0, 1.0,
                                 8, 0.2, 1.0))
    , den(std::make_unique<TH2D>("den",
                                 "Denominator; q_{inv}; k_{T} (GeV);",
                                 200, 0.0, 1.0,
                                 8, 0.2, 1.0))
    { }

  void Loop()
    {
      Event ev;

      for (;;) {
        {
          std::unique_lock<std::mutex> lk(mutex);
          cv.wait(lk, [&] { return !event_queue.empty(); });
          auto ev_ptr = std::move(event_queue.front());
          if (ev_ptr == nullptr) {
            event_queue.pop_front();
            break;
          }
          ev = std::move(*ev_ptr);
          event_queue.pop_front();
        }

        FillRealPairs(ev);
        FillMixedPairs(ev);

        mixing_queue.emplace_back(std::move(ev));
        if (mixing_queue.size() > 5) {
          mixing_queue.pop_front();
        }
      }
    }

  void FillRealPairs(const Event &ev)
    {
      for (UInt_t i=0; i<ev.particles.size(); ++i) {
        for (UInt_t j=i+1; j<ev.particles.size(); ++j) {
          const auto &pa = ev.particles[i],
                     &pb = ev.particles[j];
          FillPair(pa, pb, *num);
        }
      }
    }

  void FillMixedPairs(const Event &ev)
    {
      for (const auto &revent : mixing_queue) {
        for (const auto &pa : revent.particles) {
          for (const auto &pb : ev.particles) {
            FillPair(pa, pb, *den);
          }
        }
      }
    }

  void FillPair(const Particle &pa,
                const Particle &pb,
                TH2D &h)
    {
      auto q = pa.p - pb.p;
      auto kT = 0.5 * (pa.p.Pt() + pb.p.Pt());
      h.Fill(kT, -q.Mag());
    }

  std::thread spawn()
    {
      return std::thread(&FemtoRunner::Loop, this);
    }
};
